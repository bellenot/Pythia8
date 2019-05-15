// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

#include "Analysis.h"
#include "Basics.h"
#include "Beams.h"
#include "Event.h"
#include "HadronLevel.h"
#include "Information.h"
#include "LesHouches.h"
#include "MultipleInteractions.h"
#include "ParticleData.h"
#include "ParticleDecays.h"
#include "PartonDistributions.h"
#include "PartonLevel.h"
#include "PhaseSpace.h"
#include "ProcessLevel.h"
#include "Pythia6.h"
#include "Settings.h"
#include "SigmaProcess.h"
#include "Stdlib.h"

namespace Pythia8 {
 
//**************************************************************************

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor allows to specify default output stream. 
  Pythia(ostream& osIn = cout) : osDefault(osIn) { 
    // Initial values for parton density functions (PDF's).
    pdfAPtr = 0; pdfBPtr = 0; pdfAnew = false; pdfBnew = false;
    // Read in files with all flags, modes and parameters.
    settings.init();
    // Read in files with all particle data.
    particleData.init();
    // Write the Pythia banner to output. 
    banner();
    // Default for some flags.
    hasPythia6 = false;
    hasLHA = false;
    isInit = false;
  } 

  // Destructor to undo the PDF's created with new.
  ~Pythia() { if (pdfAnew) delete pdfAPtr; if (pdfBnew) delete pdfBPtr; } 

  // Read in one update for flag/mode/parameter/Pythia6 from a single line.
  bool readString(string, bool warn = true); 
 
  // Read in updates for flag/mode/parameter/Pythia6 from user-defined file.
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
 
  // Generate the next event.
  bool next(); 

  // Main routine to provide final statistics on generation.
  void statistics();

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

private: 

  // Static initialization data, normally only set once.
  static bool afterProcessLevel, afterPartonLevel, checkEvent;
  static int nErrList;
  static double epTolerance;

  // Constants: could only be changed in the code itself.
  static const int NTRY;

  // Default stream for output. Does not yet work??
  ostream& osDefault;

  // Write the Pythia banner, with symbol and version information.
  void banner();

  // Initialization routine to set up kinematic and more.
  bool init(bool inCMframeIn);

  // Initialization routine for all accessible static data members.
  void initStatic();

  // Check that the final event makes sense.
  bool check();

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
