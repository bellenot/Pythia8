// This file contains the main class for performing resonance decays.
// ResonanceDecays: handles the sequential decay of resonances in process.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_ResonanceDecays_H
#define Pythia8_ResonanceDecays_H

#include "Basics.h"
#include "Event.h"
#include "Information.h"
#include "ParticleData.h"
#include "ProcessContainer.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {
  
//**************************************************************************

// The ResonanceDecays class handles the sequential decay of resonances
// that are part of the hard process (t, W, Z, H, SUSY,...).

class ResonanceDecays {

public:

  // Constructor. 
  ResonanceDecays() {} 
 
  // Initialize static data members.
  // static void initStatic();
 
  // Generate the next decay sequence.
  bool next( Event& process); 

private: 

  // Pick decay channel of resonance.
  // bool pickChannel();

  // Select trial kinemantics in resonance rest frame.
  // bool pickKinematics();

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ResonanceDecays_H
