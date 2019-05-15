// This file contains the main class for parton-level event generation
// PartonLevel: administrates showers, multiple interactions and remnants.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_PartonLevel_H
#define Pythia8_PartonLevel_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "Beams.h"
#include "TimeShower.h"
#include "SpaceShower.h"
#include "MultipleInteractions.h"

namespace Pythia8 {
 
//**************************************************************************

// The PartonLevel class contains the top-level routines to generate
// the partonic activity of an event.

class PartonLevel {

public:

  // Constructor. 
  PartonLevel() {} 
 
  // Initialization assuming all necessary data already read.
  bool init( BeamParticle& beamA, BeamParticle& beamB, int strategyIn = 0);
 
  // Generate the next parton-level process.
  bool next( BeamParticle& beamA, BeamParticle& beamB, Event& process, 
    Event& event); 

  // Print statistics, if any.
  void statistics();

private: 

  // Initialization data, normally only set once.
  bool ISR, MI, FSRinProcess, FSRinResonances;

  // Constants: could only be changed in the code itself.
  static const int NTRY;

  // The generator class for timelike showers
  TimeShower times;

  // The generator class for spacelike showers
  SpaceShower space;

  // The generator class for multiple interactions.
  MultipleInteractions multi;

  // The generator class to construct beam-remnant kinematics. 
  BeamRemnants remnants;

  // Set up the hard process, excluding subsequent resonance decays.
  void setupHardSys( BeamParticle& beamA, BeamParticle& beamB, 
    Event& process, Event& event);
  // Keep track of how much of hard process has been handled.
  int nHardDone;

  // Set up the hard process, special case if all partons already given.
  bool setupSimpleSys( Event& process, Event& event);

  // Perform showers in resonance decay chains.
  void resonanceShowers( Event& process, Event& event); 
  // Position in main event record of hard partons before showers.
  vector<int> iPosBefShow;

  // Event generation strategy.
  int strategyLHA;
  
};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_PartonLevel_H
