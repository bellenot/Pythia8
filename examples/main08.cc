// File: main08.cc
// This is a simple test program. 
// It illustrates how event files can be used in Pythia8.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"
using namespace Pythia8; 
int main() {

  int nPrint = 2;                              // Number of events to print.
  Pythia pythia;                               // Generator.
  pythia.readString("PartonLevel:MI = off");   // No multiple interactions.
  pythia.readString("SpaceShower:pTmin = 1.0"); // Change pTmin cutoff of ISR.
  LHAinitPythia6 lhaInit("ttsample.init");     // Les Houches initialization object.
  LHAevntPythia6 lhaEvnt("ttsample.evnt");     // Les Houches event object.
  pythia.init(&lhaInit, &lhaEvnt);             // Initialize with pointers.
  cout << lhaInit;                             // List initialization information.
  Hist nFinal("final particle multiplicity",100,-0.5,499.5);        // Histogram.
  
  int iEvent = 0;                              // Begin event loop 
  while (pythia.next()) {                      // Generate event until none left.
    if (iEvent++ < nPrint) {                   // List first few events.
      cout << lhaEvnt;                         // List Les Houches input event.
      pythia.process.list();                   // List Pythia hard-process event.
      pythia.event.list();                     // List Pythia complete event.
    }                                          // End listing.  
    int nFin = 0;                              // Sum up final multiplicity
    for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].remains()) nFin++;
    nFinal.fill(nFin);                         // Fill histogram.
  }                                            // End of event loop.

  cout << nFinal;                              // Print histogram.
  return 0;                                    // Done.
}
