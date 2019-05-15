// File: main07.cc
// This is a simple test program. 
// It illustrates how Pythia6 processes can be used in Pythia8.
// All input is specified in the main07.cmnd file.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"

using namespace Pythia8; 

int main() {

  // Generator. Shorthand for the event and the (static) Settings.
  Pythia pythia;
  Event& event = pythia.event;
  Settings& settings = pythia.settings;

  // Read in commands from external file.
  pythia.readFile("main07.cmnd");    

  // Extract settings to be used in the main program.
  int idBeamA = settings.mode("Main:idBeamA");
  int idBeamB = settings.mode("Main:idBeamB");
  double eCM = settings.parameter("Main:eCM");
  int nEvent = settings.mode("Main:numberOfEvents");
  int nList = settings.mode("Main:numberToList");
  int nShow = settings.mode("Main:timesToShow");
  bool showChangedSettings = settings.flag("Main:showChangedSettings");
  bool showAllSettings = settings.flag("Main:showAllSettings");

  // Initialization for Pythia6 event input.
  pythia.init( idBeamA, idBeamB, eCM);

  // List changed data.
  if (showChangedSettings) settings.listChanged();
  if (showAllSettings) settings.listAll();

  // Histograms.
  double epTol = 1e-6 * eCM;
  Hist epCons("deviation from energy-momentum conservation",100,0.,epTol);
  Hist nFinal("final particle multiplicity",100,-0.5,799.5);
  Hist dnparticledy("dn/dy for particles",100,-10.,10.);

  // Begin event loop.
  int nPace = max(1,nEvent/nShow); 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%nPace == 0) cout << " Now begin event " << iEvent << "\n";

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.process.list();
      event.list();
    }

    // Loop over final particles in the event. 
    int nFin = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) if (event[i].remains()) {
      nFin++;
      pSum += event[i].p();
      dnparticledy.fill(event[i].y());
    }

    // Check and print event with too big energy-momentum deviation.
    nFinal.fill(nFin);
    double epDev = abs(pSum.e() - eCM) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);
    if (epDev > epTol) {
      cout << " Warning! Event with epDev = " << scientific 
           << setprecision(4) << epDev << " now listed:";
      event.list();
    }

  // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.statistics();
  cout << epCons << nFinal << dnparticledy; 

  return 0;
}
