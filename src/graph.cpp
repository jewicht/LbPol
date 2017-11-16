//g++ -o graph graph.cxx `root-config --libs --cflags` -lRooFit
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH3F.h"
#include "TLatex.h"
#include "TLegend.h"

#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"

#include "rootlogon.h"
#include "functions.h"
//#include "functions-roofit.cxx"
//include "progressbar.cxx"
//include "weighting.cxx"


void multiplyweights(RooDataSet* data, RooRealVar& wtot, const TString& waccname, const TString& wmassname, bool normalizetowmass, TString wcorrname)
{
  const int nEntries=data->numEntries();
  if (nEntries<=0) return;
  RooDataSet dataww6("dataww6","dataww6",RooArgList(wtot));
  
  const double defaultval=-666.;
  const double defaultwcorrval=-666.;

  double sumofwmass=0.;
  double sumofwacc=0.;
  double sumofwtot=0.;
  //  double sumofweights=0.;

  const bool debug=false;

  double wtotvec[nEntries];
  double wmassval,waccval;
  double wcorrval=1.;
  
  bool haswcorr=false;
  if (wcorrname!="") haswcorr=true;

  const RooArgSet* blu = data->get(0);
  wmassval=blu->getRealValue(wmassname, defaultval);;
  waccval=blu->getRealValue(waccname,defaultval);
  if (wmassval==defaultval) {
    std::cout << "No " << wmassname << " in dataset" << std::endl;
    exit(1);
  }
  if (waccval==defaultval) {
    std::cout << "No " << waccname << " in dataset" << std::endl;
    exit(1);
  }
  // if (haswcorr) {
  //   wcorrval=blu->getRealValue(wcorrname,defaultval);
  //   if (wcorrval==defaultval) {
  //     std::cout << "No " << wcorrname << " in dataset" << std::endl;
  //     exit(1);
  //   }
  // }
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    wmassval=blu->getRealValue(wmassname, defaultval);;
    waccval=blu->getRealValue(waccname,defaultval);
    if (haswcorr) wcorrval=blu->getRealValue(wcorrname,defaultwcorrval);
    if (debug) std::cout << "wacc = " << waccval << " ; wmass = " << wmassval << " ; wcorr = " << wcorrval <<std::endl;

    wtotvec[i]= waccval*wmassval*wcorrval;

    sumofwmass+=wmassval;
    sumofwacc+=waccval;
    sumofwtot+=wtotvec[i];
  }

  //  if (debug) mygetchar;
  //rw weights
  // make that sum (wtot) == sum(wmass)
  for (int i=0;i<nEntries;i++) { 
    if (normalizetowmass) wtotvec[i]*=sumofwmass/sumofwtot;
    wtot.setVal(wtotvec[i]);
    dataww6.add(RooArgSet(wtot));
  }

  data->merge(&dataww6);
  
}

struct myvars
{
  TString name, title;
};

RooRealVar mass("M","M(#LambdaJ/#psi)",5500.,5750.,"MeV/c^{2}");
RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);
RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);
RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);
RooRealVar wmass("wmass","wmass",0.);
RooRealVar wacc("wacc","wacc",0.);
RooRealVar wcorr("wcorr","wcorr",0.);
RooRealVar wtot("wtot","wtot",0.);
RooRealVar cat("cat","cat",0.);

int main()
{
  lhcbstyle();
  gStyle->SetPadBottomMargin(0.2);
  TCanvas c("c","c",100,100);

  std::vector<myvars> datavec;
  myvars tmp;
  
  //  tmp.name="dd/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetAll-sim-data.txt";
  //  tmp.title="";
  //  datavec.push_back(tmp);
  tmp.name="Lb-dd-down/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetDown-sim-data.txt";
  tmp.title="Lb-DD-down";
  datavec.push_back(tmp);
  tmp.name="Lb-dd-up/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetUp-sim-data.txt";
  tmp.title="Lb-DD-up";
  datavec.push_back(tmp);
  tmp.name="Lbbar-dd-down/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetDown-sim-data.txt";
  tmp.title="Lbbar-DD-down";
  datavec.push_back(tmp);
  tmp.name="Lbbar-dd-up/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetUp-sim-data.txt";
  tmp.title="Lbbar-DD-up";
  datavec.push_back(tmp);
  // tmp.name="Lbbar-ll-down/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetDown-sim-data.txt";
  // tmp.title="";
  // datavec.push_back(tmp);
  // tmp.name="Lbbar-ll-up/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetUp-sim-data.txt";
  // tmp.title="";
  // datavec.push_back(tmp);
  // tmp.name="Lb-ll-down/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetDown-sim-data.txt";
  // tmp.title="";
  // datavec.push_back(tmp);
  // tmp.name="Lb-ll-up/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetUp-sim-data.txt";
  // tmp.title="";
  // datavec.push_back(tmp);
  tmp.name="ll/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetAll-sim-data.txt";
  tmp.title="LL";
  datavec.push_back(tmp);
  // tmp.name="tot/datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetAll-sim-data.txt";
  // tmp.title="TOT";
  // datavec.push_back(tmp);


  std::vector<RooDataSet*> data;
  //  TH1D* histo[3][2][2][2];
  //  memset(data,0,8*3*sizeof(RooDataSet*));
  std::vector<TH1D*> coshist[3];


  for (unsigned int i=0;i<datavec.size();i++) {

    RooDataSet* data0 = RooDataSet::read(datavec[i].name,RooArgList(cat,mass,wmass,wacc,costheta,costheta1,costheta2));
    multiplyweights(data0,wtot,wmass.GetName(),wacc.GetName(),false,"");
    RooDataSet* data1 = new RooDataSet("data","data",data0,*data0->get(),"","wtot");
    data.push_back(data1);
    TH3F histo3d("histo3d","histo3d",
		 25,-1.,1.,
		 25,-1.,1.,
		 25,-1.,1.
		 );
    data1->fillHistogram(&histo3d,RooArgList(costheta,costheta1,costheta2));
    coshist[0].push_back(histo3d.ProjectionX("projx" + istr(i)));
    coshist[1].push_back(histo3d.ProjectionY("projy" + istr(i)));
    coshist[2].push_back(histo3d.ProjectionZ("projz" + istr(i)));
    delete data0;
  }

  //  c.SaveAs("blublu-costheta1.eps");
  //  return 1;
  std::cout << __LINE__ << std::endl;

  costheta.setBins(10);
  costheta1.setBins(10);
  costheta2.setBins(10);


  for (int k=0;k<3;k++) {
    TLegend leg(0.6,0.7,0.9,0.9);
    leg.SetFillColor(kWhite);

    for (unsigned int i=0;i<data.size();i++) {
      //      coshist[k][i]->SetLineColor( 1+ i );
      //    std::cout << histo1d[j] << std::endl;
      //    histo1d[j]->Draw();//Normalized("SAME");

      switch (i) {
      case 0: coshist[k][i]->SetLineColor(kBlue); coshist[k][i]->SetFillColor(kBlue);  break;
      case 1: coshist[k][i]->SetLineColor(kRed); coshist[k][i]->SetFillColor(kRed); coshist[k][i]->SetFillStyle(3345); break;
      case 2: coshist[k][i]->SetLineColor(kBlack); break;
      case 3: coshist[k][i]->SetLineColor(kGreen); break;
      case 4: coshist[k][i]->SetLineColor(kMagenta); break;
      }

      if (i>0) {
	if (i==1) {
	  coshist[k][i]->DrawNormalized("HIST SAME");
	} else {
	  coshist[k][i]->DrawNormalized("SAME");//Normalized("SAME");
	}
      } else {
	coshist[k][i]->SetMinimum(0.);
	coshist[k][i]->SetMaximum(coshist[k][i]->GetMaximum()*2.0);
	if (k==0) coshist[k][i]->GetXaxis()->SetTitle("cos#theta");
	if (k==1) coshist[k][i]->GetXaxis()->SetTitle("cos#theta_{1}");
	if (k==2) coshist[k][i]->GetXaxis()->SetTitle("cos#theta_{2}");
	coshist[k][i]->DrawNormalized("HIST");//Normalized();
      }
      leg.AddEntry(coshist[k][i], datavec[i].title);
    }

    leg.Draw();
    c.SaveAs("blublu-costheta" +istr(k)+".eps");

  }


  return 1;



 

  std::vector<myvars> vars;
  //  myvars tmp;
  tmp.name="ap2";  tmp.title="|a_{+}|^{2}";  vars.push_back(tmp);
  tmp.name="am2";  tmp.title="|a_{-}|^{2}";  vars.push_back(tmp);
  tmp.name="bp2";  tmp.title="|b_{+}|^{2}";  vars.push_back(tmp);
  tmp.name="bm2";  tmp.title="|b_{-}|^{2}";  vars.push_back(tmp);

  tmp.name="Pb";  tmp.title="P_{b}";  vars.push_back(tmp);
  tmp.name="ab";  tmp.title="#alpha_{b}";  vars.push_back(tmp);
  tmp.name="r0";  tmp.title="r_{0}";  vars.push_back(tmp);
  tmp.name="r1";  tmp.title="r_{1}";  vars.push_back(tmp);
  //  vars.push_back("am2");
  //  vars.push_back("bp2");
  //  vars.push_back("bm2");

  //  gStyle->SetOptFit(1111);

//  TCanvas c("c","c",150,100);

  TLatex* lat[6][3];
  
  const double xstart=0.20;
  const double xlength=0.68;
  const double y0=0.15;
  const double y1=0.09;
  const double y2=0.04;
  for (int i=0;i<2;i++) {
    for (int j=0;j<1;j++) {
      for (int k=0;k<2;k++) {
	TString lblbar("#Lambda^{0}_{b}");
	TString ddll("DD");
	TString updown("up");
	if (i==1) lblbar="#bar{#Lambda}_{b}^{0}";
	if (j==1) ddll="LL";
	if (k==1) updown="down";

	//	TString text("#splitline{" + lblbar + "-" + ddll + "}{" + updown + "}");
	TString text("#splitline{#splitline{" + lblbar + "}{" + ddll + "}}{" + updown + "}");
	lat[i*2+j*2+k][0] = new TLatex( xstart + xlength * (i*2+j*2+k) /5., y0, lblbar );
	lat[i*2+j*2+k][1] = new TLatex( xstart + xlength * (i*2+j*2+k) /5., y1, ddll );
	lat[i*2+j*2+k][2] = new TLatex( xstart + xlength * (i*2+j*2+k) /5., y2, updown );
      }
    }
  }
  lat[4][0] = new TLatex(xstart + xlength *4./5., y0, "");
  lat[4][1] = new TLatex(xstart + xlength *4./5., y1, "LL");
  lat[4][2] = new TLatex(xstart + xlength *4./5., y2, "");
  lat[5][0] = new TLatex(xstart + xlength *5./5., y0, "");
  lat[5][1] = new TLatex(xstart + xlength *5./5., y1, "ALL");
  lat[5][2] = new TLatex(xstart + xlength *5./5., y2, "");
  
  for (int j=0; j<6 ; j++ ) {
    for (int k=0; k<3; k++) {
      lat[j][k]->SetNDC();
      lat[j][k]->SetTextSize(0.05);
    }
  }
  for (unsigned int i=0; i<vars.size(); i++) {
    TGraphErrors ap2graph("graph-" + vars[i].name + ".txt","%lg %lg %lg" );

    ap2graph.Fit("pol0","","",0.5,8.5);

    ap2graph.GetXaxis()->SetTitle("");
    ap2graph.GetYaxis()->SetTitle(vars[i].title);
    ap2graph.GetXaxis()->SetLabelOffset(1000);

    ap2graph.Draw("A*");

    for (int j=0; j<6 ; j++ ) for (int k=0; k<3; k++) lat[j][k]->Draw();

    c.SaveAs("graph-" + vars[i].name + ".eps");
  }
  

}
