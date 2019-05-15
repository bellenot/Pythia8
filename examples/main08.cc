// File: main08.cc
// This is a simple test program. 
// It illustrates how to combine subruns in pT bins.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"

using namespace Pythia8; 

int main() {

  // Number of events to generate per bin, and to list.
  int nEvent = 10000;
  int nList = 1;

  // Book histograms.
  Hist pTraw("pTHat distribution, unweighted", 100, 0., 1000.);
  Hist pTnorm("pTHat distribution, weighted", 100, 0., 1000.);
  Hist pTpow4("pTHat distribution, pT4*weighted", 100, 0., 1000.);
  Hist pTpow6("pTHat distribution, pT6*weighted", 100, 0., 1000.);
  Hist pTnormPart("pTHat distribution, weighted", 100, 0., 1000.);
  Hist pTpow4Part("pTHat distribution, pT4*weighted", 100, 0., 1000.);
  Hist pTpow6Part("pTHat distribution, pT6*weighted", 100, 0., 1000.);

  // Generator.
  Pythia pythia;

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings = pythia.settings;
  Info& info = pythia.info;

  // Set up to generate QCD jets.
  pythia.readString("HardQCD:all = on");  
  pythia.readString("Pythia:afterProcessLevel = off");  

  // Set up five pT bins - last one open-ended.
  double pTlimit[6] = {100., 150., 250., 400., 600., 0.};
  for (int iBin = 0; iBin < 5; ++iBin) {
     settings.parm("PhaseSpace:pTHatMin", pTlimit[iBin]);  
     settings.parm("PhaseSpace:pTHatMax", pTlimit[iBin + 1]);  

     // Initialize for LHC.
     pythia.init( 2212, 2212, 14000.);

    // List changed data.
    settings.listChanged();

    // Reset local histograms (that need to be rescaled before added).
    pTnormPart.null();
    pTpow4Part.null();
    pTpow6Part.null();

    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Quit if failure.
      if (!pythia.next()) break;
 
      // List first few events, only hard process.
      if (iEvent < nList) pythia.process.list();

      // Fill hard scale of event.
      double pTHat = info. pTHat();
      pTraw.fill( pTHat );
      pTnormPart.fill( pTHat );
      pTpow4Part.fill( pTHat, pow4(pTHat) );
      pTpow6Part.fill( pTHat, pow3(pTHat*pTHat) );


    // End of event loop. Statistics.
    }
    pythia.statistics();

    // Normalize each case to cross section/(bin * event), and add to sum.
    double sigmaNorm = info.sigmaGen() / (10. * nEvent);
    pTnorm += sigmaNorm * pTnormPart;
    pTpow4 += sigmaNorm * pTpow4Part;
    pTpow6 += sigmaNorm * pTpow6Part;

  // End of pT-bin loop.
  }

  // Output histograms.
  cout << pTraw << pTnorm << pTpow4 << pTpow6; 

  // Done.
  return 0;
}
