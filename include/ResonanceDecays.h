// This file contains the main class for performing resonance decays.
// ResonanceDecays: handles the sequential decay of resonances in process.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_ResonanceDecays_H
#define Pythia8_ResonanceDecays_H

#include "Basics.h"
#include "Event.h"
#include "Information.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "ResonanceProperties.h"
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

  // Constants: could only be changed in the code itself.
  static const int    NTRYDECAY;
  static const double MSAFETY;

  // Select colours in decay.
  bool pickColours(Particle& decayer, Event& process);

  // Select trial kinemantics in resonance rest frame.
  void pickKinematics();

  // Flavour, colour and momentum information.
  int    id0, id0Abs, id1, id2, col0, acol0, col1, acol1, col2, acol2, 
    colType1, colType2;
  double m0, m1, m2;
  Vec4   p1, p2;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ResonanceDecays_H
