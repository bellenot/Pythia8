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
#include "SusyLesHouches.h"

namespace Pythia8 {
 
//**************************************************************************

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor. (See Pythia.cc file.)
  Pythia();

  // Destructor. (See Pythia.cc file.)
  ~Pythia();

  // Read in one update for setting/particledata/Pythia6 from a single line.
  bool readString(string, bool warn = true); 
 
  // Read in updates for setting/particledata/Pythia6 from user-defined file.
  bool readFile(string, bool warn = true);

  // Possibility to pass in pointers to PDF's. Usage optional. 
  bool PDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn);

  // Possibility to pass in pointer for external handling of some decays.
  bool decayPtr( DecayHandler* decayHandlePtrIn, 
    vector<int> handledParticles) 
    { return hadronLevel.decayPtr( decayHandlePtrIn, handledParticles);}  

  // Possibility to pass in pointer for external random number generation.
  bool rndmEnginePtr( RndmEngine* rndmPtrIn) 
    { return Rndm::rndmEnginePtr( rndmPtrIn);}   

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
  static bool partonLevelOn, hadronLevelOn, checkEvent;
  static int nErrList;
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

  // Pointers to the parton distributions of the two incoming beams.
  PDF* pdfAPtr;  
  PDF* pdfBPtr; 

  // Auxiliary to set parton densities among list of possibilities.
  PDF* setPDFPtr(int idIn);

  // Keep track when "new" has been used and needs a "delete".  
  bool pdfAnew, pdfBnew;

  // The two incoming beams.
  BeamParticle beamA;
  BeamParticle beamB;

  // Should Pythia6 be used for generating external events?
  bool hasPythia6;

  // LHAinit and LHAevnt objects for generating external events.
  bool hasLHA;
  LHAinit* lhaInitPtr;
  LHAevnt* lhaEvntPtr;
  int strategyLHA;

  // The main generator class to define the core process of the event.
  ProcessLevel processLevel;

  // The main generator class to produce the parton level of the event.
  PartonLevel partonLevel;

  // The main generator class to produce the hadron level of the event.
  HadronLevel hadronLevel;

  // ErrorMsg is a static class, so not needed here, except as reminder.
  ErrorMsg errorMsg; 

  // Properties found at the initialization of the event generator.
  bool isInit, inCMframe;
  int idA, idB;  
  double mA, mB, eA, eB, pzA, pzB, eCM, betaZ, gammaZ;

  // information for error checkout.
  int nErrEvent;
  vector<int> iErrId, iErrNan;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
