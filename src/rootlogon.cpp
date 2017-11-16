#include <iostream>

#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "rootlogon.h"



void rootlogon()
{
  std::cout << std::endl;
  std::cout << "executing rootlogon.C:" << std::endl;
  std::cout << std::endl;
  //cout << "setting Thomas' personal plot style, modified by Luc-san" << endl;
  //cout << "                      " << endl;
  
  ////////////////////////////////////////////////////////////////////
  // this style gives reasonably nice looking histograms as long
  // as they are more or less quadratic. For very long histograms,
  // adjustements are needed. For instance, for a canvas with 1x5
  // histograms:
  //  TCanvas* c1 = new TCanvas("c1", "L0 muons", 600, 800);
  //  c1->Divide(1,5);
  // adaptions like the following will be needed:
  //  gStyle->SetTickLength(0.05,"x");
  //  gStyle->SetTickLength(0.01,"y");
  //  gStyle->SetLabelSize(0.15,"x");
  //  gStyle->SetLabelSize(0.1,"y");
  //  gStyle->SetStatW(0.15);
  //  gStyle->SetStatH(0.5);
  ////////////////////////////////////////////////////////////////////
  
  TStyle *myStyle= new TStyle("myStyle","personal plots style");
  
  // font: 130 gives much nicer screen and gif look, but then 
  // I cannot have vertical plot labels (for y-axis) nor 
  // greek symbols etc. (TLatex) - therefore stick to precision 2
  // and use 132...
  //Int_t myFont = 132;
  Int_t myFont = 62;
  Int_t myFontforLegend = 42;
  
  // use plain black on white colors
  myStyle->SetFrameBorderMode(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  //myStyle->SetTitleColor(0);
  //myStyle->SetStatColor(0);
  //  myStyle->SetFillColor(19); !! laisser en commentaire
  myStyle->SetPalette(1);
  
  
  // set the paper & margin sizes
  myStyle->SetPaperSize(20,20);
  myStyle->SetPadTopMargin(0.05);
  myStyle->SetPadRightMargin(0.05); // increase for colz plots!!
  myStyle->SetPadBottomMargin(0.16);
  //myStyle->SetPadLeftMargin(0.12);
  myStyle->SetPadLeftMargin(0.20);
  
  // use large Times-Roman fonts
  myStyle->SetTextFont(myFontforLegend);
  //myStyle->SetTextSize(0.08);
  myStyle->SetTextSize(0.10);
  myStyle->SetLabelFont(myFont,"x");
  myStyle->SetLabelFont(myFont,"y");
  myStyle->SetLabelFont(myFont,"z");
  //myStyle->SetLabelSize(0.05,"x");
  //myStyle->SetLabelSize(0.05,"y");
  //myStyle->SetLabelSize(0.05,"z");
  myStyle->SetLabelSize(0.055,"x");
  myStyle->SetLabelSize(0.055,"y");
  myStyle->SetLabelSize(0.055,"z");
  myStyle->SetTitleFont(myFont);
  //myStyle->SetTitleSize(0.06,"x");
  //myStyle->SetTitleSize(0.06,"y");
  //myStyle->SetTitleSize(0.06,"z");
  myStyle->SetTitleSize(0.07,"x");
  myStyle->SetTitleSize(0.07,"y");
  myStyle->SetTitleSize(0.07,"z");

  myStyle->SetTitleOffset(1.4, "y");
  
  // use bold lines and markers

  //mod feb 10 2012
  //  myStyle->SetMarkerStyle(8);
  myStyle->SetMarkerSize(0.1);//myStyle->SetMarkerSize(0.5);
  myStyle->SetHistLineWidth(4);//was 2
  myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  
  // by default, only display statistics box:
  myStyle->SetOptStat(1110);  // show only nent, mean, rms
  myStyle->SetOptStat(0000);  // show only nent, mean, rms
  myStyle->SetOptTitle(0);
  myStyle->SetOptFit(0);
  
  // look of the statistics box:
  myStyle->SetStatBorderSize(1);
  myStyle->SetStatFont(myFont);
  myStyle->SetStatFontSize(0.05);
  myStyle->SetStatX(0.95);
  myStyle->SetStatY(0.95);
  
  // put tick marks on top and RHS of plots
  //myStyle->SetPadTickX(1);
  //myStyle->SetPadTickY(1);
  myStyle->SetPadTickX(0);
  myStyle->SetPadTickY(0);
  
  // myStyle->SetHatchesSpacing(1.);
  // myStyle->SetHatchesLineWidth(2);
  
  
  myStyle->SetLineWidth(4);
  
  gROOT->SetStyle("myStyle");
  gROOT->ForceStyle();
  
}


// all users - please change the name of this file to lhcbStyle.C
// Commits to lhcbdocs svn of .C files are not allowed

// Rev 136730

void lhcbstyle()
{

  // // define names for colours
  // Int_t black  = 1;
  // Int_t red    = 2;
  // Int_t green  = 3;
  // Int_t blue   = 4;
  // Int_t yellow = 5; 
  // Int_t magenta= 6;
  // Int_t cyan   = 7;
  // Int_t purple = 9;
  

////////////////////////////////////////////////////////////////////
// PURPOSE:
//
// This macro defines a standard style for (black-and-white) 
// "publication quality" LHCb ROOT plots. 
//
// USAGE:
//
// Include the lines
//   gROOT->ProcessLine(".L lhcbstyle.C");
//   lhcbStyle();
// at the beginning of your root macro.
//
// Example usage is given in myPlot.C
//
// COMMENTS:
//
// Font:
// 
// The font is chosen to be 132, this is Times New Roman (like the text of
//  your document) with precision 2.
//
// "Landscape histograms":
//
// The style here is designed for more or less square plots.
// For longer histograms, or canvas with many pads, adjustements are needed. 
// For instance, for a canvas with 1x5 histograms:
//  TCanvas* c1 = new TCanvas("c1", "L0 muons", 600, 800);
//  c1->Divide(1,5);
//  Adaptions like the following will be needed:
//  gStyle->SetTickLength(0.05,"x");
//  gStyle->SetTickLength(0.01,"y");
//  gStyle->SetLabelSize(0.15,"x");
//  gStyle->SetLabelSize(0.1,"y");
//  gStyle->SetStatW(0.15);
//  gStyle->SetStatH(0.5);
//
// Authors: Thomas Schietinger, Andrew Powell, Chris Parkes, Niels Tuning
// Maintained by Editorial board member (currently Niels)
///////////////////////////////////////////////////////////////////

  // // Use times new roman, precision 2 
  // Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // lhcbFont = 42;
  // // Line thickness
  // Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // // Text size
  // Double_t lhcbTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
  
  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
  //  lhcbStyle->SetPadBottomMargin(0.16);
  //  lhcbStyle->SetPadLeftMargin(0.14);
  lhcbStyle->SetPadBottomMargin(0.20);
  lhcbStyle->SetPadLeftMargin(0.20);
  
  //for acceptance plots
  //  lhcbStyle->SetPadLeftMargin(0.25);

  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  //  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  //  lhcbStyle->SetFrameLineWidth(1.);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  // // JW  lhcbStyle->SetMarkerSize(1.0);
  lhcbStyle->SetMarkerSize(0.1);

  //JW
  lhcbStyle->SetEndErrorSize(0.);

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  //  lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0); 
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
  
  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  // // add LHCb label
  // TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
  //                          0.87 - gStyle->GetPadTopMargin(),
  //                          gStyle->GetPadLeftMargin() + 0.20,
  //                          0.95 - gStyle->GetPadTopMargin(),
  //                          "BRNDC");
  // lhcbName->AddText("LHCb");
  // lhcbName->SetFillColor(0);
  // lhcbName->SetTextAlign(12);
  // lhcbName->SetBorderSize(0);

  // TText *lhcbLabel = new TText();
  // lhcbLabel->SetTextFont(lhcbFont);
  // lhcbLabel->SetTextColor(1);
  // lhcbLabel->SetTextSize(lhcbTSize);
  // lhcbLabel->SetTextAlign(12);

  // TLatex *lhcbLatex = new TLatex();
  // lhcbLatex->SetTextFont(lhcbFont);
  // lhcbLatex->SetTextColor(1);
  // lhcbLatex->SetTextSize(lhcbTSize);
  // lhcbLatex->SetTextAlign(12);

  std::cout << "-------------------------" << std::endl;  
  std::cout << "Set LHCb Style - Feb 2012" << std::endl;
  std::cout << "-------------------------" << std::endl;  
  
}
