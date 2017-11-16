# This is the Makefile to make the macros 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

ROOFITLIBS   := -lRooFit -lRooFitCore -lRooStats -lMinuit2 -lHistPainter

CXX           = g++
#CXXFLAGS      = -g -Wall -I..
#LDFLAGS      = -g -Wall
CXXFLAGS      =  -fPIC -O2 -Wall -Iinclude/
LDFLAGS      = -O2 -Wall
LD            = g++

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = -Lrootmod/lib $(ROOTLIBS) $(SYSLIBS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

# Keep this defined if you measure efficiencies down to 15 GeV. Comment it
# out if you are measuring efficiences down to 25. (It should match what
# EMCert_x was run with.) It basically controls some defaults in the fitting
# code for the invariant mass shape. The name comes because em_cert originally
# used 25 GeV as the lower cutoff, but when top-style cuts were implemented,
# the pT cut was dropped to 15 GeV.
#CXXFLAGS     += -DDO_TOP_CUTS

#*****************************************************************************
#* Targets
#*****************************************************************************

OBJECT1 = rootlogon.o computew3.o computew5.o RooLbtoJpsiL0PDF5.o RooLbtoJpsiL0PDF5-dict.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o RooLbtoJpsiL0wAccPDF3.o RooLbtoJpsiL0wAccPDF3-dict.o    RooLegendre5Pdfv2.o RooLegendre5Pdfv2-dict.o createRooLegendre3.o createRooLegendre5.o RooLegendre3Pdfv2.o RooLegendre3Pdfv2-dict.o RooChebychev3Pdf.o RooChebychev3Pdf-dict.o calcacceptanceclass.o findminmaxtf3.o functions-roofit.o progressbar.o test.o 
BINARY1 = test

OBJECT2 = computew5.o computew3.o writedecayfile.o
BINARY2 = writedecayfile

OBJECT3 = makepdf.o 
BINARY3 = makepdf

OBJECT4 = splittree.o
BINARY4 = splittree

OBJECT5 = rootlogon.o  RooLegendre5Pdfv2.o RooLegendre5Pdfv2-dict.o RooLegendre3Pdfv2.o RooLegendre3Pdfv2-dict.o RooChebychev3Pdf.o RooChebychev3Pdf-dict.o computew3.o  RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o RooLbtoJpsiL0wAccPDF3.o RooLbtoJpsiL0wAccPDF3-dict.o createRooLegendre3.o createRooLegendre5.o findminmaxtf3.o  functions-roofit.o progressbar.o testlegendrepdf.o
BINARY5 = testlegendrepdf

OBJECT6 = rootlogon.o   functions-roofit.o  progressbar.o cuts.o compmcdata.o
BINARY6 = compmcdata

OBJECT7 = calcangles.o progressbar.o angles.o
BINARY7 = calcangles

OBJECT8 = reducedata.o cuts.o progressbar.o
BINARY8 = reducedata

OBJECT9 = rootlogon.o weighting.o cuts.o functions-roofit.o progressbar.o computew5.o RooLbtoJpsiL0PDF5.o RooLbtoJpsiL0PDF5-dict.o computew3.o computew3v2.o computew3v3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o RooLbtoJpsiL0PDF3v2.o RooLbtoJpsiL0PDF3v2-dict.o RooLbtoJpsiL0PDF3v3.o RooLbtoJpsiL0PDF3v3-dict.o RooLbtoJpsiL0PDF3v4.o RooLbtoJpsiL0PDF3v4-dict.o calcacceptanceclass.o ConfigFile.o  RooLbtoJpsiL0wAccPDF3.o RooLbtoJpsiL0wAccPDF3-dict.o createRooLegendre3.o createRooLegendre2.o RooChebychev3Pdf.o RooChebychev3Pdf-dict.o RooLegendre3Pdfv2.o RooLegendre3Pdfv2-dict.o RooLegendre2Pdf.o RooLegendre2Pdf-dict.o  angles.o datafit.o
BINARY9 = datafit

#OBJECT10 = findminmaxtf3.o   calcacceptanceclass.o calcacceptance.o 
#BINARY10 = calcacceptance

OBJECT11 = trainbdt.o cuts.o
BINARY11 = trainbdt

OBJECT12 = calcbdt.o progressbar.o cuts.o
BINARY12 = calcbdt

OBJECT14 = optbdt.o functions-roofit.o progressbar.o cuts.o rootlogon.o 
BINARY14 = optbdt

OBJECT15 = calcksmass.o
BINARY15 = calcksmass

OBJECT16 = cuts.o calcacceptanceclass.o rootlogon.o  RooLegendre3Pdfv2.o RooLegendre3Pdfv2-dict.o RooLegendre2Pdf.o RooLegendre2Pdf-dict.o RooChebychev3Pdf.o RooChebychev3Pdf-dict.o computew3.o RooLbtoJpsiL0wAccPDF3.o RooLbtoJpsiL0wAccPDF3-dict.o createRooLegendre3.o findminmaxtf3.o  functions-roofit.o progressbar.o weighting.o ConfigFile.o fitacceptance3.o  createRooLegendre2.o
BINARY16 = fitacceptance3

OBJECT17 = convolution.o rootlogon.o computew3.o  RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o 
BINARY17 = convolution

OBJECT18 = cuts.o rootlogon.o progressbar.o functions-roofit.o reweightmc.o
BINARY18 = reweightmc

OBJECT19 = savetruemc.o
BINARY19 = savetruemc

OBJECT20 = addmissingvardata.o
BINARY20 = addmissingvardata

OBJECT21 = rootlogon.o  calcacceptanceclass.o computew3v3.o RooLbtoJpsiL0PDF3v3.o RooLbtoJpsiL0PDF3v3-dict.o  computew3v2.o RooLbtoJpsiL0PDF3v2.o RooLbtoJpsiL0PDF3v2-dict.o computew3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o functions-roofit.o progressbar.o ConfigFile.o weighting.o toy.o cuts.o
BINARY21 = toy

OBJECT22 = rootlogon.o  functions-roofit.o progressbar.o RooLegendre5Pdfv2.o RooLegendre5Pdfv2-dict.o  createRooLegendre5.o findminmaxtf3.o cuts.o fitacceptance5.o calcacceptanceclass5.o
BINARY22 = fitacceptance5

OBJECT23 = rootlogon.o  functions-roofit.o progressbar.o cuts.o fitmcmass.o 
BINARY23 = fitmcmass

OBJECT24 = computew3.o writedecayfile3.o
BINARY24 = writedecayfile3

OBJECT25 = rootlogon.o  calcacceptanceclass.o computew3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o functions-roofit.o progressbar.o ConfigFile.o weighting.o createRooLegendre3.o RooLegendre3Pdfv2-dict.o  RooLegendre3Pdfv2.o toyprod.o RooLbtoJpsiL0wAccPDF3.o RooLbtoJpsiL0wAccPDF3-dict.o RooChebychev3Pdf.o RooChebychev3Pdf-dict.o
BINARY25 = toyprod

OBJECT26 = rootlogon.o functions-roofit.o progressbar.o computew3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o compgaussroofit.o
BINARY26 = compgaussroofit

OBJECT27 = rootlogon.o functions-roofit.o progressbar.o computew3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o fitevtgen.o
BINARY27 = fitevtgen

OBJECT28 = graph.o rootlogon.o 
BINARY28 = graph

OBJECT29 = fitbias.o functions-roofit.o  progressbar.o rootlogon.o weighting.o calcacceptanceclass.o
BINARY29 = fitbias

OBJECT30 = angularres.o functions-roofit.o  progressbar.o rootlogon.o computew3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o cuts.o
BINARY30 = angularres

OBJECT31 = fitaccphi1phi2.o functions-roofit.o progressbar.o rootlogon.o RooLegendre2Pdf.o RooLegendre2Pdf-dict.o createRooLegendre2.o cuts.o ConfigFile.o
BINARY31 = fitaccphi1phi2

OBJECT32 = convandintegrate.o  rootlogon.o
BINARY32 = convandintegrate

OBJECT33 = rootlogon.o weighting.o cuts.o functions-roofit.o progressbar.o computew5.o RooLbtoJpsiL0PDF5.o RooLbtoJpsiL0PDF5-dict.o computew3.o computew3v2.o computew3v3.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o RooLbtoJpsiL0PDF3v2.o RooLbtoJpsiL0PDF3v2-dict.o RooLbtoJpsiL0PDF3v3.o RooLbtoJpsiL0PDF3v3-dict.o RooLbtoJpsiL0PDF3v4.o RooLbtoJpsiL0PDF3v4-dict.o calcacceptanceclass.o ConfigFile.o  RooLbtoJpsiL0wAccPDF3.o RooLbtoJpsiL0wAccPDF3-dict.o createRooLegendre3.o createRooLegendre2.o RooChebychev3Pdf.o RooChebychev3Pdf-dict.o RooLegendre3Pdfv2.o RooLegendre3Pdfv2-dict.o RooLegendre2Pdf.o RooLegendre2Pdf-dict.o  angles.o mctest.o
BINARY33 = mctest

OBJECT34 = mergegraph.o  rootlogon.o
BINARY34 = mergegraph

OBJECT35 = plotpdf.o  rootlogon.o RooLbtoJpsiL0PDF3.o RooLbtoJpsiL0PDF3-dict.o computew3.o  functions-roofit.o progressbar.o
BINARY35 = plotpdf


######## OLD UNUPDATED SOURCE FILES  ###############
# (removed from standard buid; add them to the 
#  BINARIES definition below if you want to use them)
####### all binaries #############

BINARIES = $(BINARY2)  $(BINARY3)  $(BINARY4)   $(BINARY6)   $(BINARY7)   $(BINARY8) ${BINARY9}  ${BINARY10} ${BINARY11} ${BINARY12} ${BINARY14} ${BINARY15} ${BINARY16} ${BINARY17} ${BINARY18} ${BINARY19} ${BINARY20} ${BINARY21} ${BINARY22} ${BINARY23} ${BINARY24} ${BINARY25} ${BINARY26} ${BINARY27} ${BINARY28} ${BINARY29} ${BINARY30} ${BINARY31} ${BINARY32} ${BINARY33} ${BINARY34} ${BINARY35}


#*****************************************************************************
#* Dependencies
#*****************************************************************************

all: $(BINARIES)

$(BINARY1) :   $(OBJECT1)
	$(LD) $(LDFLAGS) $(OBJECT1) $(LIBS) $(ROOFITLIBS) -o $(BINARY1)

$(BINARY2) :   $(OBJECT2)
	$(LD) $(LDFLAGS) $(OBJECT2) $(LIBS) -o $(BINARY2)

$(BINARY3) :   $(OBJECT3)
	$(LD) $(LDFLAGS) $(OBJECT3) $(LIBS) $(ROOFITLIBS) -o $(BINARY3)

$(BINARY4) :   $(OBJECT4)
	$(LD) $(LDFLAGS) $(OBJECT4) $(LIBS) -o $(BINARY4)

$(BINARY5) :   $(OBJECT5)
	$(LD) $(LDFLAGS) $(OBJECT5) $(LIBS) $(ROOFITLIBS) -o $(BINARY5)

$(BINARY6) :   $(OBJECT6)
	$(LD) $(LDFLAGS) $(OBJECT6) $(LIBS) $(ROOFITLIBS) -o $(BINARY6)

$(BINARY7) :   $(OBJECT7)
	$(LD) $(LDFLAGS) $(OBJECT7) $(LIBS) -o $(BINARY7)

$(BINARY8) :   $(OBJECT8)
	$(LD) $(LDFLAGS) $(OBJECT8) $(LIBS) -o $(BINARY8)

$(BINARY9) :   $(OBJECT9)
	$(LD) $(LDFLAGS) $(OBJECT9) $(LIBS)  $(ROOFITLIBS) -o $(BINARY9)

$(BINARY10) :   $(OBJECT10)
	$(LD) $(LDFLAGS) $(OBJECT10) $(LIBS) -o $(BINARY10)

$(BINARY11) :   $(OBJECT11)
	$(LD) $(LDFLAGS) $(OBJECT11) $(LIBS) -lTMVA -o $(BINARY11)

$(BINARY12) :   $(OBJECT12)
	$(LD) $(LDFLAGS) $(OBJECT12) $(LIBS) $(ROOFITLIBS) -lTMVA -o $(BINARY12)

$(BINARY13) :   $(OBJECT13)
	$(LD) $(LDFLAGS) $(OBJECT13) $(LIBS) -lGui -o $(BINARY13)

$(BINARY14) :   $(OBJECT14)
	$(LD) $(LDFLAGS) $(OBJECT14) $(LIBS)  $(ROOFITLIBS) -o $(BINARY14)

$(BINARY15) :   $(OBJECT15)
	$(LD) $(LDFLAGS) $(OBJECT15) $(LIBS) -o $(BINARY15)

$(BINARY16) :   $(OBJECT16)
	$(LD) $(LDFLAGS) $(OBJECT16) $(LIBS)  $(ROOFITLIBS) -o $(BINARY16)

$(BINARY17) :   $(OBJECT17)
	$(LD) $(LDFLAGS) $(OBJECT17) $(LIBS)  $(ROOFITLIBS) -o $(BINARY17)

$(BINARY18) :   $(OBJECT18)
	$(LD) $(LDFLAGS) $(OBJECT18) $(LIBS) $(ROOFITLIBS)  -o $(BINARY18)

$(BINARY19) :   $(OBJECT19)
	$(LD) $(LDFLAGS) $(OBJECT19) $(LIBS) -o $(BINARY19)

$(BINARY20) :   $(OBJECT20)
	$(LD) $(LDFLAGS) $(OBJECT20) $(LIBS) -o $(BINARY20)

$(BINARY21) :   $(OBJECT21)
	$(LD) $(LDFLAGS) $(OBJECT21) $(LIBS)  $(ROOFITLIBS) -o $(BINARY21)

$(BINARY22) :   $(OBJECT22)
	$(LD) $(LDFLAGS) $(OBJECT22) $(LIBS)  $(ROOFITLIBS) -o $(BINARY22)

$(BINARY23) :   $(OBJECT23)
	$(LD) $(LDFLAGS) $(OBJECT23) $(LIBS)  $(ROOFITLIBS) -o $(BINARY23)

$(BINARY24) :   $(OBJECT24)
	$(LD) $(LDFLAGS) $(OBJECT24) $(LIBS) -o $(BINARY24)

$(BINARY25) :   $(OBJECT25)
	$(LD) $(LDFLAGS) $(OBJECT25) $(LIBS)  $(ROOFITLIBS) -o $(BINARY25)

$(BINARY26) :   $(OBJECT26)
	$(LD) $(LDFLAGS) $(OBJECT26) $(LIBS)  $(ROOFITLIBS) -o $(BINARY26)

$(BINARY27) :   $(OBJECT27)
	$(LD) $(LDFLAGS) $(OBJECT27) $(LIBS)  $(ROOFITLIBS) -o $(BINARY27)

$(BINARY28) :   $(OBJECT28)
	$(LD) $(LDFLAGS) $(OBJECT28) $(LIBS) $(ROOFITLIBS) -o $(BINARY28)

$(BINARY29) :   $(OBJECT29)
	$(LD) $(LDFLAGS) $(OBJECT29) $(LIBS) $(ROOFITLIBS) -o $(BINARY29)

$(BINARY30) :   $(OBJECT30)
	$(LD) $(LDFLAGS) $(OBJECT30) $(LIBS) $(ROOFITLIBS) -o $(BINARY30)

$(BINARY31) :   $(OBJECT31)
	$(LD) $(LDFLAGS) $(OBJECT31) $(LIBS) $(ROOFITLIBS) -o $(BINARY31)

$(BINARY32) :   $(OBJECT32)
	$(LD) $(LDFLAGS) $(OBJECT32) $(LIBS) -o $(BINARY32)

$(BINARY33) :   $(OBJECT33)
	$(LD) $(LDFLAGS) $(OBJECT33) $(LIBS) $(ROOFITLIBS) -o $(BINARY33)

$(BINARY34) :   $(OBJECT34)
	$(LD) $(LDFLAGS) $(OBJECT34) $(LIBS) -o $(BINARY34)

$(BINARY35) :   $(OBJECT35)
	$(LD) $(LDFLAGS) $(OBJECT35) $(LIBS) $(ROOFITLIBS) -o $(BINARY35)


RooLbtoJpsiL0PDF3-dict.cpp: include/RooLbtoJpsiL0PDF3.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLbtoJpsiL0PDF3v2-dict.cpp: include/RooLbtoJpsiL0PDF3v2.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLbtoJpsiL0PDF3v3-dict.cpp: include/RooLbtoJpsiL0PDF3v3.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLbtoJpsiL0PDF3v4-dict.cpp: include/RooLbtoJpsiL0PDF3v4.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLbtoJpsiL0wAccPDF3-dict.cpp: include/RooLbtoJpsiL0wAccPDF3.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^	

RooLbtoJpsiL0PDF5-dict.cpp: include/RooLbtoJpsiL0PDF5.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

#RooLegendre5Pdf-dict.cpp: include/RooLegendre5Pdf.h
#	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLegendre5Pdfv2-dict.cpp: include/RooLegendre5Pdfv2.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

# RooLegendre3Pdf-dict.cpp: include/RooLegendre3Pdf.h
# 	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLegendre3Pdfv2-dict.cpp: include/RooLegendre3Pdfv2.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooLegendre2Pdf-dict.cpp: include/RooLegendre2Pdf.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

RooChebychev3Pdf-dict.cpp: include/RooChebychev3Pdf.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

Tetris-dict.cpp: Tetris.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^


%.o: src/%.cpp 
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
%.o: roofitmodels/%.cpp 
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
.PHONY : clean

clean:
	$(RM) *.o *-dict.* *.pcm $(BINARIES)


