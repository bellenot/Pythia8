// This file contains the main class for process-level event generation.
// ProcessLevel: administrates the selection of "hard" process.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_ProcessLevel_H
#define Pythia8_ProcessLevel_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "Beams.h"
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
  bool init( int idA, int idB, double eCM, bool hasLHAin = false, 
    LHAinit* lhaInitPtrIn = 0, LHAevnt* lhaEvntPtrIn = 0);
 
  // Generate the next "hard" process.
  bool next( BeamParticle& beamA, BeamParticle& beamB, Event& process); 

  // Main routine to provide final statistics on generation.
  // void statistics();

private: 

  // LHAinit and LHAevnt objects for generating external events.
  bool hasLHA;
  LHAinit* lhaInitPtr;
  LHAevnt* lhaEvntPtr;
  int strategyLHA;

  // Initialize event generation from Pythia 6.3.
  void initPythia6( int idA, int idB, double eCM);

  // Read in the hard process from the Les Houches Accord.
  bool getLHAevnt( BeamParticle& beamA, BeamParticle& beamB, Event& process);

  // Read in the hard process, special case if all partons already given.
  bool getSimpleLHAevnt( Event& process);

  // Add any junctions to the process event record list.
  void findJunctions( Event& process);

  // Reordered list of LHA partons.
  vector<int> newPos;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ProcessLevel_H
