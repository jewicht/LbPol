#include <math.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"

#include "rootlogon.h"

int main(int argc, char** argv)
{
  lhcbstyle();
  const int vecsize=50;
  
  double vecx[vecsize];
  double vecy[vecsize];
  for (int j=0; j<vecsize; j++) {vecx[j]=0.;vecy[j]=0.;}

  for (int i=2;i<argc;i++) {
    TFile inputfile(argv[i]);
    TGraph* graph = (TGraph*)inputfile.Get("Graph");
    
    for (int j=0; j<vecsize; j++) {
      double x,y;
      graph->GetPoint(j,x,y);
if (isnan(y) or (y!=y) or isinf(y)) {
        std::cerr << "Problem in " << argv[i] << std::endl;
        //return 1;
      }
      vecx[j]=x;
      vecy[j]+=y;//exp(-y);
    }
    inputfile.Close();
  }
  for (int j=0; j<vecsize; j++) {
    vecy[j]/=(argc-1);
//    vecy[j]=-log(vecy[j]);
  }


  TGraph outputgraph(vecsize,vecx,vecy);
  outputgraph.SetName("Graph");

  TFile outputfile(argv[1], "recreate");
  outputgraph.Write();
  outputfile.Close();

  // TCanvas c("c","c",100,100);
  // outputgraph.Draw("AC*");
  // c.SaveAs("outputgraph.eps");

}
