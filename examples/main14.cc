// File: main14.cc
// This is a simple test program. 
// It illustrates how to use internal Pythia8 processes,
// with special emphasis on hard processes.
// All input is specified in the main14.cmnd file.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"

using namespace Pythia8; 

//**************************************************************************

int main() {

  // Generator. Shorthand for the event and for settings.
  Pythia pythia;
  Event& event = pythia.event;
  Settings& settings = pythia.settings;

  // Read in commands from external file.
  pythia.readFile("main14.cmnd");    

  // Extract settings to be used in the main program.
  int idBeamA = settings.mode("Main:idBeamA");
  int idBeamB = settings.mode("Main:idBeamB");
  bool inCMframe = settings.flag("Main:inCMframe");
  double eCM = settings.parameter("Main:eCM");
  double eBeamA = settings.parameter("Main:eBeamA");
  double eBeamB = settings.parameter("Main:eBeamB");
  int nEvent = settings.mode("Main:numberOfEvents");
  int nList = settings.mode("Main:numberToList");
  int nShow = settings.mode("Main:timesToShow");
  int nAbort = settings.mode("Main:timesAllowErrors");
  bool showChangedSettings = settings.flag("Main:showChangedSettings");
  bool showAllSettings = settings.flag("Main:showAllSettings");
  bool showChangedParticleData 
    = settings.flag("Main:showChangedParticleData");
  bool showAllParticleData = settings.flag("Main:showAllParticleData");
 
  // Initialization for Pythia6 event input.
  if (inCMframe) pythia.init( idBeamA, idBeamB, eCM);
  else pythia.init( idBeamA, idBeamB, eBeamA, eBeamB);

  // List changed data.
  if (showChangedSettings) settings.listChanged();
  if (showAllSettings) settings.listAll();

  // List particle data.  
  if (showChangedParticleData) ParticleDataTable::listChanged();
  if (showAllParticleData) ParticleDataTable::listAll();

  // Book histograms.
  Hist pThard("scattering pT spectrum", 100, 0., 100.); 
  Hist nInt("number of interactions", 100, -0.5, 99.5);
  Hist enhanceFac("b enhancement factor", 100, 0., 10.);
  Hist nChg("number of charged particles", 100, -0.5, 199.5);
  Hist pTchg("charged particle pT spectrum", 100, 0., 10.); 
 
  // Begin event loop.
  int nShowPace = max(1,nEvent/nShow); 
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%nShowPace == 0) cout << " Now begin event " 
      << iEvent << endl;

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.info.list();
      pythia.process.list();
      event.list();
    }

    // Hard process transverse momentum.
    double pT = pythia.process[5].pT();
    // Alternative expression below only works for internal processes. 
    // double pT = pythia.info.pTHat();
    pThard.fill( pT );

    // Number of interactions and enhancement from impact parameter.
    int nMI = pythia.info.nMI();
    nInt.fill( nMI );
    double enhance = pythia.info.enhanceMI();
    enhanceFac.fill( enhance );

    // Charged multiplicity and pT spectrum.
    int nch = 0;
    for (int i = 1; i < event.size(); ++i)
      if (event[i].remains() && event[i].isCharged()) {
        ++nch; 
        pTchg.fill( event[i].pT() );
      }
    nChg.fill( nch );

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.statistics();
  cout << pThard << nInt << enhanceFac << nChg << pTchg; 

  // Done.
  return 0;
}
