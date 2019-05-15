// This file contains the main class for process-level event generation.
// ProcessLevel: administrates the selection of "hard" process.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_ProcessLevel_H
#define Pythia8_ProcessLevel_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "Beams.h"
#include "Event.h"
#include "Information.h"
#include "ProcessContainer.h"
#include "LesHouches.h"
#include "Pythia6.h"

namespace Pythia8 {
  
//**************************************************************************

// The ProcessLevel class contains the top-level routines to generate
// the characteristic "hard" process of an event.

class ProcessLevel {

public:

  // Constructor. 
  ProcessLevel() {} 
 
  // Initialization assuming all necessary data already read.
  bool init( Info& info, PDF* pdfAPtr = 0, PDF* pdfBPtr = 0, 
    bool hasPythia6In = false, bool hasLHAin = false, 
    LHAinit* lhaInitPtrIn = 0, LHAevnt* lhaEvntPtrIn = 0);
 
  // Generate the next "hard" process.
  bool next( Info& info, Event& process); 

  // Print statistics on cross sections and number of events.
  void statistics(ostream& = cout);

private: 

  // Which machinery is used to generate events?
  bool hasInternal, hasPythia6, hasLHA;

  // Vector of containers of internally-generated processes.
  vector<ProcessContainer*> containerPtrs;
  double sigmaMaxSum;

  // LHAinit and LHAevnt objects for generating external events.
  LHAinit* lhaInitPtr;
  LHAevnt* lhaEvntPtr;
  int strategyLHA;

  // SigmaTotal object needed to handle soft QCD processes.
  SigmaTotal sigmaTot;

  // Initialize the internal event generation machinery.
  void initInternal( Info info, PDF* pdfAPtr, PDF* pdfBPtr,
    ostream& = cout);

  // Initialize event generation from Pythia 6.3.
  void initPythia6( int idA, int idB, double eCM);

  // Generate the next internal event.
  bool getInternalEvnt( Info& info, Event& process);

  // Read in the hard process from the Les Houches Accord.
  bool getLHAevnt( Info& info, Event& process);

  // Read in the hard process, special case if all partons already given.
  bool getSimpleLHAevnt( Event& process);

  // Add any junctions to the process event record list.
  void findJunctions( Event& process);

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ProcessLevel_H
