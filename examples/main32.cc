// main32.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Mikhail Kirsanov, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how HepMC can be interfaced to Pythia8.
// It takes input from main32.cmnd and puts events on a HepMC file,
// with all analysis intended to happen afterwards.

#include "Pythia.h"

#include "HepMCInterface.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_Ascii.h"
//#include "HepMC/IO_AsciiParticles.h"

using namespace Pythia8; 
int main() {

  // Interface for conversion from Pythia8::Event to HepMC one. 
  HepMC::I_Pythia8 ToHepMC;
  //  ToHepMC.set_crash_on_problem();

  // Specify file where HepMC events will be stored.
  HepMC::IO_Ascii ascii_io("hepmcout32.dat",std::ios::out);
  // HepMC::IO_AsciiParticles ascii_io("hepmcout32.dat",std::ios::out);

  // Generator. 
  Pythia pythia;

  // Read in commands from external file.
  pythia.readFile("main32.cmnd");    

  // Extract settings to be used in the main program.
  int    idBeamA   = pythia.mode("Main:idBeamA");
  int    idBeamB   = pythia.mode("Main:idBeamB");
  double eCM       = pythia.parm("Main:eCM");
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nShow     = pythia.mode("Main:timesToShow");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");
  bool   showCS    = pythia.flag("Main:showChangedSettings");
  bool   showAS    = pythia.flag("Main:showAllSettings");
  bool   showCPD   = pythia.flag("Main:showChangedParticleData");
  bool   showAPD   = pythia.flag("Main:showAllParticleData");
 
  // Initialization.
  pythia.init( idBeamA, idBeamB, eCM);

  // List settings.
  if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  // Begin event loop.
  int nShowPace = max(1,nEvent/nShow); 
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%nShowPace == 0) cout << " Now begin event " 
      << iEvent << endl;

    // Generate event. Skip if erroneous. Quit if too many failures.   
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    // Convert event record to HepMC format and output to file.
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia.event, hepmcevt );
    ascii_io << hepmcevt;
    delete hepmcevt;

  // End of event loop. Statistics. 
  }
  pythia.statistics();

  // Done.
  return 0;
}
