// File: main12.cc
// This is a simple test program. 
// It compares Pythia6 and Pythia8 cross sections for many processes.
// All input is specified in the main12.cmnd file, where the real action is.
// Copyright C 2007 Torbjorn Sjostrand

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
  double eCM = settings.parm("Main:eCM");
  double eBeamA = settings.parm("Main:eBeamA");
  double eBeamB = settings.parm("Main:eBeamB");
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
  Hist mHard("scattering mass spectrum", 100, 0., 1000.); 
  Hist pThard("scattering pT spectrum", 100, 0., 100.); 
  Hist mHigh("max mass 3,4", 100, 0., 200.); 
  Hist mLow("min mass 3,4", 100, 0., 200.); 
 
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

    // Hard process mass.
    double m = pythia.info.mHat();
    mHard.fill(m);

    // Hard process transverse momentum - only meaningful for 2 -> 2.
    if (pythia.info.nFinal() == 2) { 
      double pT = pythia.info.pTHat();
      pThard.fill( pT );

      // Product masses - only meaningful for 2 -> 2.
      double m3Hat = pythia.process[5].m();
      double m4Hat = pythia.process[6].m();
      mHigh.fill(max(m3Hat,m4Hat));
      mLow.fill(min(m3Hat,m4Hat));
    }

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.statistics();
  cout << mHard << pThard << mHigh << mLow ;

  // Done.
  return 0;
}
