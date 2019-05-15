// File: main02.cc
// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.
// Copyright C 2007 Torbjorn Sjostrand
#include "Pythia.h"
using namespace Pythia8; 
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");    
  pythia.readString("PhaseSpace:mHatMin = 80.");    
  pythia.readString("PhaseSpace:mHatMax = 120.");    
  pythia.init( 2212, -2212, 1960.);
  Hist pTZ("dN/dpTZ", 100, 0., 100.);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) {pythia.info.list(); pythia.event.list();}
    // Loop over particles in event. Find last Z0 copy. Fill its pT. 
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].id() == 23) iZ = i;
    pTZ.fill( pythia.event[iZ].pT() );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.statistics();
  cout << pTZ; 
  return 0;
}
