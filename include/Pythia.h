// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

#include "Analysis.h"
#include "Basics.h"
#include "Event.h"
#include "HadronLevel.h"
#include "Information.h"
#include "LesHouches.h"
#include "PartonLevel.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "ProcessLevel.h"
#include "PythiaStdlib.h"
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
  bool readFile(string, bool warn = true);

  // Possibility to pass in pointers to PDF's. Usage optional. 
  bool setPDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn, PDF* pdfHardAPtrIn = 0, 
    PDF* pdfHardBPtrIn = 0);

  // Possibility to pass in pointer for external handling of some decays.
  bool setDecayPtr( DecayHandler* decayHandlePtrIn, 
    vector<int> handledParticles) 
    { return hadronLevel.decayPtr( decayHandlePtrIn, handledParticles);}  

  // Possibility to pass in pointer for external random number generation.
  bool setRndmEnginePtr( RndmEngine* rndmEnginePtrIn) 
    { return Rndm::rndmEnginePtr( rndmEnginePtrIn);}  

  // Possibility to pass in pointer for user hooks. Usage optional. 
  bool setUserHooksPtr( UserHooks* userHooksPtrIn) 
    { userHooksPtr = userHooksPtrIn; return true;} 

  // Possibility to pass in pointer for external showers. Usage optional. 
  bool setShowerPtr( TimeShower* timesDecPtrIn, 
    TimeShower* timesPtrIn = 0, SpaceShower* spacePtrIn = 0) 
    { timesDecPtr = timesDecPtrIn; timesPtr = timesPtrIn;
    spacePtr = spacePtrIn; return true;} 

  // Initialization with two beams specified.
  bool init( int idAin, int idBin, double eAin, double eBin);

  // Initialization in the CM frame.
  bool init( int idAin, int idBin, double eCMin);

  // Initialization according to the Les Houches Accord.
  bool init( LHAinit* lhaInitPtrIn, LHAevnt* lhaEvntPtrIn);

  // Initialization by a Les Houches Event File.
  bool init( string LesHouchesEventFile);

  // Initialization of data only - not enough to generate events.
  bool initData() {initStatic(); particleData.initBWmass(); return true;}
 
  // Generate the next event.
  bool next(); 

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

  // Static initialization data, normally only set once.
  static bool   doPartonLevel, doHadronLevel, checkEvent;
  static int    nErrList;
  static double epTolerance;

  // Constants: could only be changed in the code itself.
  static const int NTRY;

  // Write the Pythia banner, with symbol and version information.
  void banner(ostream& os = cout);

  // Initialization routine to set up kinematic and more.
  bool init();

  // Initialization routine for all accessible static data members.
  void initStatic();

  // Initialization routine for SUSY spectra.
  bool initSLHA();

  // Check that the final event makes sense.
  bool check(ostream& os = cout);

  // Keep track when "new" has been used and needs a "delete".  
  bool pdfAnew, pdfBnew, pdfHardNew;

  // Pointers to the parton distributions of the two incoming beams.
  PDF* pdfAPtr;  
  PDF* pdfBPtr; 

  // Extra PDF pointers to be used in hard processes only. 
  PDF* pdfHardAPtr;  
  PDF* pdfHardBPtr; 

  // Auxiliary to set parton densities among list of possibilities.
  PDF* getPDFPtr(int idIn, int sequence = 1);

  // The two incoming beams.
  BeamParticle beamA;
  BeamParticle beamB;

  // LHAinit and LHAevnt objects for generating external events.
  bool hasLHA;
  LHAinit* lhaInitPtr;
  LHAevnt* lhaEvntPtr;
  int strategyLHA;

  // Pointer to userHooks object for user interaction with program.
  UserHooks* userHooksPtr;
  bool hasUserHooks, doVetoProcess, doVetoPartons;

  // Keep track when "new" has been used and needs a "delete".  
  bool timesNew, spaceNew;

  // Pointers to timelike and spacelike showers.
  TimeShower*  timesDecPtr;
  TimeShower*  timesPtr;
  SpaceShower* spacePtr;

  // The main generator class to define the core process of the event.
  ProcessLevel processLevel;

  // The main generator class to produce the parton level of the event.
  PartonLevel partonLevel;

  // The main generator class to produce the hadron level of the event.
  HadronLevel hadronLevel;

  // ErrorMsg is a static class, so not needed here, except as reminder.
  ErrorMsg errorMsg; 

  // Properties found at the initialization of the event generator.
  bool   isConstructed, isInit, inCMframe;
  int    idA, idB;  
  double mA, mB, eA, eB, pzA, pzB, eCM, betaZ, gammaZ;

  // information for error checkout.
  int    nErrEvent;
  vector<int> iErrId, iErrNan;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
