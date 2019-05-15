// File: main12.cc
// This is a simple test program. 
// It illustrates how to use internal Pythia8 processes.
// All input is specified in the main12.cmnd file.
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
  pythia.readFile("main12.cmnd");    

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
  Hist pTspec("scattering pT spectrum", 100, 0., 100.); 
  Hist pTspecLow("scattering pT spectrum", 100, 0., 1.);
  Hist tSpecEl("elastic |t| spectrum", 100, 0., 5.);
  Hist tSpecSD("single diffractive |t| spectrum", 100, 0., 5.); 
  Hist tSpecDD("double diffractive |t| spectrum", 100, 0., 5.); 
  Hist mSpec("scattering mass spectrum", 100, 0., 100.); 
  Hist mLogSpec("log10(scattering mass spectrum)", 100, 0., 4.); 

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

    // Study event.
    double pT = pythia.process[5].pT(); 
    // double pT = pythia.info.pTHat(); // Only works for Pythia8 processes
    pTspec.fill( pT );
    pTspecLow.fill( pT );
    mSpec.fill( pythia.info.m3Hat() );
    mSpec.fill( pythia.info.m4Hat() );
    mLogSpec.fill( log10(pythia.info.m3Hat()) );
    mLogSpec.fill( log10(pythia.info.m4Hat()) );
    int code = pythia.info.code();
    double tAbs = abs(pythia.info.tHat());
    if (code == 102) tSpecEl.fill(tAbs);
    else if (code == 103 || code == 104) tSpecSD.fill(tAbs);
    else if (code == 105) tSpecDD.fill(tAbs);

  // End of event loop.
  }

  // Final statistics.
  pythia.statistics();
  cout << pTspec << pTspecLow;
  // For elastic/diffractive processes:
  // cout << tSpecEl << tSpecSD << tSpecDD << mSpec << mLogSpec;

  // Done.
  return 0;
}
