// This file contains the main class for process-level event generation.
// ProcessLevel: administrates the selection of "hard" process.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_ProcessLevel_H
#define Pythia8_ProcessLevel_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Information.h"
#include "LesHouches.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "ProcessContainer.h"
#include "PythiaStdlib.h"
#include "ResonanceDecays.h"
#include "Settings.h"
#include "UserHooks.h"

namespace Pythia8 {
  
//**************************************************************************

// The ProcessLevel class contains the top-level routines to generate
// the characteristic "hard" process of an event.

class ProcessLevel {

public:

  // Constructor. 
  ProcessLevel() {} 
 
  // Initialization assuming all necessary data already read.
  bool init( Info* infoPtrIn, BeamParticle* beamAPtrIn, 
    BeamParticle* beamBPtrIn, bool doLHAin = false, 
    LHAinit* lhaInitPtrIn = 0, LHAevnt* lhaEvntPtrIn = 0, 
    UserHooks* userHooksPtrIn = 0);
 
  // Generate the next "hard" process.
  bool next( Event& process); 

  // Accumulate and update statistics (after possible user veto).
  void accumulate();

  // Print statistics on cross sections and number of events.
  void statistics(ostream& os = cout);

private: 

  // Generic info for process generation.
  bool   doInternal, doLHA, doSecondHard, allHardSame, noneHardSame, 
         someHardSame;
  int    nImpact, startColTag2;
  double sigmaND, sumImpactFac, sum2ImpactFac;

  // Vector of containers of internally-generated processes.
  vector<ProcessContainer*> containerPtrs;
  int    iContainer;
  double sigmaMaxSum;

  // Ditto for optional choice of a second hard process.
  vector<ProcessContainer*> container2Ptrs;
  int    i2Container;
  double sigma2MaxSum;

  // Pointer to various information on the generation.
  Info* infoPtr;

  // Pointers to the two incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Pointers to LHAinit and LHAevnt for generating external events.
  LHAinit* lhaInitPtr;
  LHAevnt* lhaEvntPtr;
  int strategyLHA;

  // Pointer to userHooks object for user interaction with program.
  UserHooks* userHooksPtr;

  // SigmaTotal object needed to handle soft QCD processes.
  SigmaTotal sigmaTot;

  // ResonanceDecay object does sequential resonance decays.
  ResonanceDecays resonanceDecays;

  // Initialize information on resonances.
  bool initResonances();

  // Initialize the internal event generation machinery.
  bool initInternal( ostream& os = cout);

  // Generate the next internal event with one interaction.
  bool nextInternal( Event& process);

  // Generate the next internal event with two hard interactions.
  bool next2Internal( Event& process);

  // Append the second to the first process list.
  void combineProcessRecords( Event& process,  Event& process2);

  // Read in the hard process from the Les Houches Accord.
  bool nextLHA( Event& process);

  // Read in the hard process, special case if all partons already given.
  bool nextSimpleLHA( Event& process);

  // Add any junctions to the process event record list.
  void findJunctions( Event& process);

  // Check that colours match up.
  bool checkColours( Event& process);

  // Print statistics when two hard processes allowed.
  void statistics2(ostream& os = cout);

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ProcessLevel_H
