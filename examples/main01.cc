// File: main01.cc
// This is a simple test program. It fits on one slide in a talk. 
// It studies the charged multiplicity distribution at the LHC.
// Copyright C 2006 Torbjorn Sjostrand
#include "Pythia.h"
using namespace Pythia8; 
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("HardQCD:all = on");    
  pythia.readString("PhaseSpace:pTHatMin = 20.");    
  pythia.init( 2212, 2212, 14000.);
  Hist mult("charged multiplicity", 100, -0.5, 799.5);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) {pythia.info.list(); pythia.event.list();} 
    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) 
        ++nCharged; 
    mult.fill( nCharged );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.statistics();
  cout << mult; 
  return 0;
}
