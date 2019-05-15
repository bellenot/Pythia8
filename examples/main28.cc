// main28.cc is a part of the PYTHIA event generator.
// Copyright (C) 2009 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a test program for the extra dimensions processes.
// Author: Stefan Ask (Stefan DOT Ask AT cern DOT ch)
// Documentation: S. Ask, arXiv:0809.4750

#include "Pythia.h"

using namespace Pythia8; 

//+++ The main program.
int main() {

  //+++ Test cases
  //+++ 1 = Jet + G (real G emission) 
  //+++ 2 = Jet + U (real U emission) 
  //+++ 3 = Z + G   (real G emission) 
  //+++ 4 = Z + U   (real U emission)
  //+++ 5 = gamma gamma (LED G* exchange)
  //+++ 6 = l lbar      (LED G* exchange). 
  //+++     Note: charged leptons only!
  //+++ 7 = G* (RS resonance)
  int n_testproc = 1;  

  //+++ Number of events to generate and to list. Max number of errors.
  int nEvent     = 100;      
  int nList      = 1;       
  int nAbort     = 50;         

  //+++ Pythia generator.
  Pythia pythia;

  //+++ PYTHIA paramters:
  pythia.readString("PhaseSpace:showViolation = off");

  //+++ Test case parameters
  if (n_testproc == 1) { 
    pythia.readString("ExtraDimensionsLED:monojet = on");
    pythia.readString("ExtraDimensionsLED:n = 4");
    pythia.readString("ExtraDimensionsLED:MD = 4000.");
    pythia.readString("ExtraDimensionsLED:CutOffmode = 3");
    pythia.readString("ExtraDimensionsLED:t = 2");
    pythia.readString("5000039:m0 = 2500.");
    pythia.readString("5000039:mWidth = 1500.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 700.");
  } else if (n_testproc == 2){ 
    pythia.readString("ExtraDimensionsUnpart:gg2Ug = off");
    pythia.readString("ExtraDimensionsUnpart:qg2Uq = on");
    pythia.readString("ExtraDimensionsUnpart:qqbar2Ug = on");
    pythia.readString("ExtraDimensionsUnpart:spinU = 1");
    pythia.readString("ExtraDimensionsUnpart:dU = 1.2");
    pythia.readString("ExtraDimensionsUnpart:LambdaU = 1000");
    pythia.readString("ExtraDimensionsUnpart:lambda = 1.0");
    pythia.readString("ExtraDimensionsUnpart:CutOffmode = 0");
    pythia.readString("5000039:m0 = 300.");
    pythia.readString("5000039:mWidth = 500.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 700.");
  } else if (n_testproc == 3){
    pythia.readString("ExtraDimensionsLED:ffbar2GZ = on");
    pythia.readString("ExtraDimensionsLED:n = 6");
    pythia.readString("ExtraDimensionsLED:MD = 2000.");
    pythia.readString("ExtraDimensionsLED:CutOffmode = 1");
    pythia.readString("5000039:m0 = 4000.");
    pythia.readString("5000039:mWidth = 2500.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 50.");
  } else if (n_testproc == 4){ 
    pythia.readString("ExtraDimensionsUnpart:ffbar2UZ = on");
    pythia.readString("ExtraDimensionsUnpart:spinU = 1");
    pythia.readString("ExtraDimensionsUnpart:dU = 2.0");
    pythia.readString("ExtraDimensionsUnpart:LambdaU = 1000");
    pythia.readString("ExtraDimensionsUnpart:lambda = 1.000");
    pythia.readString("ExtraDimensionsUnpart:CutOffmode = 0");
    pythia.readString("5000039:m0 = 500.");
    pythia.readString("5000039:mWidth = 1000.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 50.");
  } else if (n_testproc == 5){ 
    pythia.readString("ExtraDimensionsLED:ffbar2gammagamma = on");
    pythia.readString("ExtraDimensionsLED:gg2gammagamma = on");
    pythia.readString("ExtraDimensionsLED:LambdaT = 3300.");
    pythia.readString("PhaseSpace:mHatMin = 800.");
  } else if (n_testproc == 6){ 
    pythia.readString("ExtraDimensionsUnpart:ffbar2llbar = on");
    pythia.readString("ExtraDimensionsUnpart:gg2llbar = off");
    pythia.readString("ExtraDimensionsUnpart:spinU = 1");
    pythia.readString("ExtraDimensionsUnpart:dU = 1.3");
    pythia.readString("ExtraDimensionsUnpart:LambdaU = 1000");
    pythia.readString("ExtraDimensionsUnpart:lambda = 1.0");
    pythia.readString("ExtraDimensionsUnpart:gXX = 0");
    pythia.readString("ExtraDimensionsUnpart:gXY = 0");
    pythia.readString("PhaseSpace:mHatMin = 300.");
  } else if (n_testproc == 7){
    pythia.readString("ExtraDimensionsG*:all = on");
  }

  //+++ Initialization for LHC.
  pythia.init( 2212, 2212, 14000.);

  //+++ Validation histograms
  Hist ETjet("dN/dETjet", 100, 0., 7000.);
  Hist mG("dN/mG", 100, 0., 7000.);

  //+++ Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%(max(1,nEvent/20)) == 0) std::cout << " Now begin event " 
						  << iEvent << "\n";
    
    //+++ Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      std::cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    //+++ List first few events. 
    if (iEvent < nList) { 
      pythia.info.list();
      pythia.process.list();
    }

    //+++ Checked particle index
    int tmp_monojet = -1;

    //+++ Particle loop
    for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {

      //+++ From hard process
      if ( pythia.event[iPart].statusAbs()  == 23 ) {

	//+++ Find graviton/unparticle
	if(pythia.event[iPart].idAbs() == 5000039){
	  mG.fill( pythia.event[iPart].m() );
	}

	//+++ Find mono-jets
	if (n_testproc == 1 || n_testproc == 2) {
	  if ( pythia.event[iPart].idAbs() <= 6 
            || pythia.event[iPart].idAbs() == 21 ){
	    if (tmp_monojet >= 0) {
	      std::cout << "More than one (hard process) mono-jet ! \n";
	    } else {
	      tmp_monojet  = iPart;
	    }
	  }
	}

      }
    }

    //+++ Validation mono-jet wrt G.Giudice et al. paper [hep-ph/9811291v2]	
    if (tmp_monojet >= 0) {
      double tmp_eta = pythia.event[tmp_monojet].eta();
      double tmp_et = pythia.event[tmp_monojet].eT();
      double tmp_et_cut = 1000;
      if ( tmp_et >=  tmp_et_cut && abs(tmp_eta) < 3 ) {
	ETjet.fill( fabs(tmp_et) );
      }    
    }
    
  }  
 
  //+++ Final statistics.
  pythia.statistics(); 

  std::cout << "-------------- Graviton mass  ------------" << "\n";
  std::cout << mG;
  std::cout << "-------------- Mono-jet check ------------" << "\n";
  std::cout << ETjet;
  std::cout << "------------------------------------------" << "\n";

  return 0;
}
