// main31.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Mikhail Kirsanov, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how HepMC can be interfaced to Pythia8.
// It studies the charged multiplicity distribution at the LHC.
// HepMC events are output to the hepmcout31.dat file.
// Written by Mikhail Kirsanov based on main01.cc.

#include "Pythia.h"

#include "HepMCInterface.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_Ascii.h"
//#include "HepMC/IO_AsciiParticles.h"

using namespace Pythia8; 
int main() {

  HepMC::I_Pythia8 ToHepMC;
  //  ToHepMC.set_crash_on_problem();

  // Specify file where HepMC events will be stored.
  HepMC::IO_Ascii ascii_io("hepmcout31.dat",std::ios::out);
  // HepMC::IO_AsciiParticles ascii_io("hepmcout31.dat",std::ios::out);

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  Event& event = pythia.event;
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

    // Convert event record to HepMC format and output to file.
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( event, hepmcevt );
    ascii_io << hepmcevt;
    delete hepmcevt;

  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.statistics();
  cout << mult; 
  return 0;
}
