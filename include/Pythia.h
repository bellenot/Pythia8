// Pythia.h is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

#include "Analysis.h"
#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "HadronLevel.h"
#include "Information.h"
#include "LesHouches.h"
#include "PartonLevel.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "ProcessLevel.h"
#include "PythiaStdlib.h"
#include "ResonanceWidths.h"
#include "Settings.h"
#include "SpaceShower.h"
#include "SusyLesHouches.h"
#include "TimeShower.h"
#include "UserHooks.h"

namespace Pythia8 {
 
//**************************************************************************

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor. (See Pythia.cc file.)
  Pythia(string xmlDir = "../xmldoc");

  // Destructor. (See Pythia.cc file.)
  ~Pythia();

  // Read in one update for a setting or particle data from a single line.
  bool readString(string, bool warn = true); 
 
  // Read in updates for settings or particle data from user-defined file.
  bool readFile(string, bool warn = true, int subrun = SUBRUNDEFAULT);
  bool readFile(string fileName, int subrun) {
    return readFile(fileName, true, subrun);}

  // Possibility to pass in pointers to PDF's.
  bool setPDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn, PDF* pdfHardAPtrIn = 0, 
    PDF* pdfHardBPtrIn = 0);

  // Possibility to pass in pointer for external handling of some decays.
  bool setDecayPtr( DecayHandler* decayHandlePtrIn, 
    vector<int> handledParticlesIn) {decayHandlePtr = decayHandlePtrIn; 
    handledParticles.resize(0); 
    for(int i = 0; i < int(handledParticlesIn.size()); ++i)
    handledParticles.push_back( handledParticlesIn[i] ); return true;}  

  // Possibility to pass in pointer for external random number generation.
  bool setRndmEnginePtr( RndmEngine* rndmEnginePtrIn) 
    { return Rndm::rndmEnginePtr( rndmEnginePtrIn);}  

  // Possibility to pass in pointer for user hooks. 
  bool setUserHooksPtr( UserHooks* userHooksPtrIn) 
    { userHooksPtr = userHooksPtrIn; return true;} 

  // Possibility to pass in pointer(s) for external cross section.
  bool setSigmaPtr( SigmaProcess* sigmaPtrIn) 
    { sigmaPtrs.push_back( sigmaPtrIn); return true;} 

  // Possibility to pass in pointer(s) for external resonance.
  bool setResonancePtr( ResonanceWidths* resonancePtrIn) 
    { resonancePtrs.push_back( resonancePtrIn); return true;} 

  // Possibility to pass in pointer for external showers.
  bool setShowerPtr( TimeShower* timesDecPtrIn, 
    TimeShower* timesPtrIn = 0, SpaceShower* spacePtrIn = 0) 
    { timesDecPtr = timesDecPtrIn; timesPtr = timesPtrIn;
    spacePtr = spacePtrIn; return true;} 

  // Initialization with two beams specified.
  bool init( int idAin, int idBin, double eAin, double eBin);

  // Initialization in the CM frame.
  bool init( int idAin, int idBin, double eCMin);

  // Initialization using the Main beam variables.
  bool init();

  // Initialization according to the Les Houches Accord.
  bool init( LHAinit* lhaInitPtrIn, LHAevnt* lhaEvntPtrIn);

  // Initialization by a Les Houches Event File.
  bool init( string LesHouchesEventFile, bool skipInit = false);
 
  // Generate the next event.
  bool next(); 

  // Special routine to allow more decays if on/off switches changed.
  bool moreDecays() {return hadronLevel.moreDecays(event);}

  // List the current Les Houches event.
  void LHAevntList(ostream& os = cout) {lhaEvntPtr->list(os);}

  // Main routine to provide final statistics on generation.
  void statistics(bool all = false);

  // Read in settings values: shorthand, not new functionality.
  bool   flag(string key) {return settings.flag(key);}
  int    mode(string key) {return settings.mode(key);} 
  double parm(string key) {return settings.parm(key);}
  string word(string key) {return settings.word(key);}

  // The event record for the parton-level central process.
  Event process;

  // The event record for the complete event history.
  Event event;

  // Information on the generation, especially current subprocess.
  Info info;

  // Settings - is static but declared here for ease of use.
  Settings settings;

  // ParticleDataTable - is static but declared here for ease of use.
  ParticleDataTable particleData;

  // SusyLesHouches - SLHA object for interface to SUSY spectra.
  SusyLesHouches slha;

private: 

  // Constants: could only be changed in the code itself.
  static const int NTRY, SUBRUNDEFAULT;

  // Initialization data, extracted from database.
  bool   doProcessLevel, doPartonLevel, doHadronLevel, checkEvent;
  int    nErrList;
  double epTolErr, epTolWarn;

  // Initialization data, extracted from init(...) call.
  bool   isConstructed, isInit, inCMframe;
  int    idA, idB;  
  double mA, mB, eA, eB, pzA, pzB, eCM, betaZ, gammaZ;

  // information for error checkout.
  int    nInitCalls, nErrEvent;
  vector<int> iErrId, iErrNan;

  // Pointers to the parton distributions of the two incoming beams.
  PDF* pdfAPtr;  
  PDF* pdfBPtr; 

  // Extra PDF pointers to be used in hard processes only. 
  PDF* pdfHardAPtr;  
  PDF* pdfHardBPtr; 

  // Keep track when "new" has been used and needs a "delete" for PDF's.  
  bool useNewPdfA, useNewPdfB, useNewPdfHard;

  // The two incoming beams.
  BeamParticle beamA;
  BeamParticle beamB;

  // LHAinit and LHAevnt objects for generating external events.
  bool doLHA, useNewLHA;
  LHAinit* lhaInitPtr;
  LHAevnt* lhaEvntPtr;

  // Pointer to external decay handler and list of particles it handles.
  DecayHandler* decayHandlePtr;
  vector<int> handledParticles;

  // Pointer to userHooks object for user interaction with program.
  UserHooks* userHooksPtr;
  bool hasUserHooks, doVetoProcess, doVetoPartons;

  // Pointers to external processes derived from the Pythia base classes.
  vector<SigmaProcess*> sigmaPtrs;  

  // Pointers to external calculation of resonance widths.
  vector<ResonanceWidths*> resonancePtrs;

  // Pointers to timelike and spacelike showers.
  TimeShower*  timesDecPtr;
  TimeShower*  timesPtr;
  SpaceShower* spacePtr;

  // Keep track when "new" has been used and needs a "delete" for showers.  
  bool useNewTimes, useNewSpace;

  // The main generator class to define the core process of the event.
  ProcessLevel processLevel;

  // The main generator class to produce the parton level of the event.
  PartonLevel partonLevel;

  // The main generator class to produce the hadron level of the event.
  HadronLevel hadronLevel;

  // ErrorMsg is a static class, so not needed here, except as reminder.
  ErrorMsg errorMsg; 

  // Write the Pythia banner, with symbol and version information.
  void banner(ostream& os = cout);

  // Check for lines in file that mark the beginning of new subrun.
  int readSubrun(string line, bool warn = true, ostream& os = cout);

  // Initialization routine to set up kinematic and more.
  bool initInternal();

  // Initialize tunes to e+e- and pp/ppbar data.
  void initTunes();

  // Initialization routine for SUSY spectra.
  bool initSLHA();

  // Add any junctions to the event record list when no ProcessLevel.
  void findJunctions();

  // Check that the final event makes sense.
  bool check(ostream& os = cout);

  // Auxiliary to set parton densities among list of possibilities.
  PDF* getPDFPtr(int idIn, int sequence = 1);

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
