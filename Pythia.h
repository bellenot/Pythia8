// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "Beams.h"
#include "ProcessLevel.h"
#include "PartonLevel.h"
#include "HadronLevel.h"
#include "LesHouches.h"
#include "Pythia6.h"
#include "ParticleDecays.h"

namespace Pythia8 {
 
//**************************************************************************

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor allows to specify default output stream. 
  Pythia(ostream& osIn = cout) : osDefault(osIn) { 
    // Initial values for parton density functions (PDF's).
    pdfAptr = 0; pdfBptr = 0; pdfAnew = false; pdfBnew = false;
    // Read in files with all flags, modes and parameters.
    settings.init();
    // Read in files with all particle data.
    particleData.init();
    // Write the Pythia banner to output. 
    banner();
    // Default for some flags.
    hasLHA = false;
  } 

  // Destructor to undo the PDF's created with new.
  ~Pythia() { if (pdfAnew) delete pdfAptr; if (pdfBnew) delete pdfBptr; } 

  // Read in one update for flag/mode/parameter/Pythia6 from a single line.
  bool readString(string, bool warn = true); 
 
  // Read in updates for flag/mode/parameter/Pythia6 from user-defined file.
  bool readFile(string, bool warn = true);

  // Possibility to pass in pointers to PDF's. Usage optional. 
  bool PDFptr( PDF* pdfAptrIn, PDF* pdfBptrIn);

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

  // Initialization for collider: pp, pbarp, ppbar, e+e-, e-e+.
  bool init( string machineIn, double eCMin);

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

  // Settings - is static but declared here for ease of use.
  Settings settings;

  // ParticleDataTable - is static but declared here for ease of use.
  ParticleDataTable particleData;

private: 

  // Static initialization data, normally only set once.
  static bool checkEvent;
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
  PDF* pdfAptr;  
  PDF* pdfBptr; 

  // Auxiliary to set parton densities among list of possibilities.
  PDF* setPDF(int idIn);

  // Keep track when "new" has been used and needs a "delete".  
  bool pdfAnew, pdfBnew;

  // The two incoming beams.
  BeamParticle beamA;
  BeamParticle beamB;

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
  bool inCMframe;
  int idA, idB;  
  double mA, mB, eA, eB, pzA, pzB, eCM, betaZ, gammaZ;

  // information for error checkout.
  int nErrEvent;
  vector<int> iErrId, iErrNan;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
