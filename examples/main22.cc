// main22.cc is a part of the PYTHIA event generator.
// Copyright (C) 2008 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates how to run SUSY processes in Pythia8.
// All input is specified in the main22.cmnd file.

#include "Pythia.h"

using namespace Pythia8; 

int main() {

  // Generator. Shorthand for the event and the (static) Settings.
  Pythia pythia;
  Event& event = pythia.event;
  Settings& settings = pythia.settings;

  // Read in commands from external file.
  pythia.readFile("main22.cmnd");    

  // Extract settings to be used in the main program.
  int nEvent = settings.mode("Main:numberOfEvents");
  int nList  = settings.mode("Main:numberToList");
  int nShow  = settings.mode("Main:timesToShow");
  bool showChangedSettings = settings.flag("Main:showChangedSettings");
  bool showAllSettings = settings.flag("Main:showAllSettings");
  bool showChangedParticleData 
    = settings.flag("Main:showChangedParticleData");
  bool showAllParticleData = settings.flag("Main:showAllParticleData");
  double eCM = settings.parm("Beams:eCM");

  // Initialize. Beam parameters set in .cmnd file.
  pythia.init();

  // List changed data.
  if (showChangedSettings) settings.listChanged();
  if (showAllSettings) settings.listAll();

  // List particle data.  
  if (showChangedParticleData) ParticleDataTable::listChanged();
  if (showAllParticleData) ParticleDataTable::listAll();

  // Histograms.
  double epTol = 1e-6 * eCM;
  Hist epCons("deviation from energy-momentum conservation",100,0.,epTol);
  Hist nFinal("final particle multiplicity",100,-0.5,799.5);
  Hist dnparticledy("dn/dy for particles",100,-10.,10.);

  // Begin event loop.
  int nPace = max(1,nEvent/nShow); 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%nPace == 0) cout << " Now begin event " << iEvent << endl;

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
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
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

/*

Pythia8060 with main16.spc (heavy squarks, msq ~ 5e9)

 | q qbar -> ~chi_10 ~chi_10                1001 |        3975        357 |   1.890e-13  4.864e-15 |
 | q qbar -> ~chi_10 ~chi_20                1002 |        6629        525 |   2.827e-13  5.940e-15 |
 | q qbar -> ~chi_10 ~chi_30                1003 |       35950       3633 |   1.855e-12  1.576e-14 |
 | q qbar -> ~chi_10 ~chi_40                1004 |        3114        210 |   1.264e-13  4.117e-15 |
 | q qbar -> ~chi_20 ~chi_20                1005 |        3162        257 |   1.390e-13  4.147e-15 |
 | q qbar -> ~chi_20 ~chi_30                1006 |      110289      10862 |   5.731e-12  2.743e-14 |
 | q qbar -> ~chi_20 ~chi_40                1007 |        2184        170 |   9.396e-14  3.580e-15 |
 | q qbar -> ~chi_30 ~chi_30                1008 |          47          5 |   1.190e-15  2.815e-16 |
 | q qbar -> ~chi_30 ~chi_40                1009 |      346486      33969 |   1.769e-11  4.793e-14 |
 | q qbar -> ~chi_40 ~chi_40                1010 |         158         12 |   7.890e-15  9.946e-16 |
 
Pythia 6 with same spectrum and neutralinos forced stable. 

 I 216 f + fbar -> ~chi1 + ~chi1    I          369          1716 I  1.822E-13 I
 I 217 f + fbar -> ~chi2 + ~chi2    I          255          1211 I  1.399E-13 I
 I 218 f + fbar -> ~chi3 + ~chi3    I            6            20 I  2.475E-15 I
 I 219 f + fbar -> ~chi4 + ~chi4    I           11            59 I  6.165E-15 I
 I 220 f + fbar -> ~chi1 + ~chi2    I          572          2824 I  2.819E-13 I
 I 221 f + fbar -> ~chi1 + ~chi3    I         3451         14329 I  1.803E-12 I
 I 222 f + fbar -> ~chi1 + ~chi4    I          233          1182 I  1.368E-13 I
 I 223 f + fbar -> ~chi2 + ~chi3    I        10808         45052 I  5.690E-12 I
 I 224 f + fbar -> ~chi2 + ~chi4    I          189           906 I  9.563E-14 I
 I 225 f + fbar -> ~chi3 + ~chi4    I        34106        142123 I  1.755E-11 I

*/
